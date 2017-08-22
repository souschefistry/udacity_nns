# https://github.com/tensorflow/tensorflow/issues/2169
# ThomasWollmann's solution of max un pooling 
# I've adapted chahld's version to handle unknown input tensor shape.

def unpool(pool, ind, ksize=[1, 2, 2, 1], scope='unpool'):
    """
       Unpooling layer after max_pool_with_argmax.
       Args:
           pool:   max pooled output tensor
           ind:      argmax indices
           ksize:     ksize is the same as for the pool
       Return:
           unpool:    unpooling tensor
    """
    with tf.variable_scope(scope):
        input_shape =  tf.shape(pool)
        output_shape = [input_shape[0], input_shape[1] * ksize[1], input_shape[2] * ksize[2], input_shape[3]]

        flat_input_size = tf.cumprod(input_shape)[-1]
        flat_output_shape = tf.stack([output_shape[0], output_shape[1] * output_shape[2] * output_shape[3]])

        pool_ = tf.reshape(pool, tf.stack([flat_input_size]))
        batch_range = tf.reshape(tf.range(tf.cast(output_shape[0], tf.int64), dtype=ind.dtype), 
                                          shape=tf.stack([input_shape[0], 1, 1, 1]))
        b = tf.ones_like(ind) * batch_range
        b = tf.reshape(b, tf.stack([flat_input_size, 1]))
        ind_ = tf.reshape(ind, tf.stack([flat_input_size, 1]))
        ind_ = tf.concat([b, ind_], 1)

        ret = tf.scatter_nd(ind_, pool_, shape=tf.cast(flat_output_shape, tf.int64))
        ret = tf.reshape(ret, tf.stack(output_shape))
        return ret
        
# test max unpooling solution from ThomasWollmann
x = tf.placeholder(tf.float32, shape=(1, None, None, None))

pool_1, argmax_1 = tf.nn.max_pool_with_argmax(x, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding='SAME', name='pool1')

pool_2, argmax_2 = tf.nn.max_pool_with_argmax(pool_1, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding='SAME', name='pool1')

pool_3, argmax_3 = tf.nn.max_pool_with_argmax(pool_2, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding='SAME', name='pool1')

pool_4, argmax_4 = tf.nn.max_pool_with_argmax(pool_3, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding='SAME', name='pool1')

pool_5, argmax_5 = tf.nn.max_pool_with_argmax(pool_4, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding='SAME', name='pool1')

# unpool(pool, ind, ksize=[1, 2, 2, 1], scope='unpool')

unpool_5 = unpool(pool_5, argmax_5, tf.shape(pool_4))

unpool_4 = unpool(unpool_5, argmax_4, tf.shape(pool_3))

unpool_3 = unpool(unpool_4, argmax_3, tf.shape(pool_2))

unpool_2 = unpool(unpool_3, argmax_2, tf.shape(pool_1))

unpool_1 = unpool(unpool_2, argmax_1, tf.shape(x))

session = tf.InteractiveSession()
session.run(tf.global_variables_initializer())
image = np.float32(cv2.imread('../data/dressthuy.JPG'))

tests = [
    pool_1.eval(feed_dict={x: [image]}),
    pool_2.eval(feed_dict={x: [image]}),
    pool_3.eval(feed_dict={x: [image]}),
    pool_4.eval(feed_dict={x: [image]}),
    pool_5.eval(feed_dict={x: [image]}),
    unpool_5.eval(feed_dict={x: [image]}),
    unpool_4.eval(feed_dict={x: [image]}),
    unpool_3.eval(feed_dict={x: [image]}),
    unpool_2.eval(feed_dict={x: [image]}),
    unpool_1.eval(feed_dict={x: [image]})
]

session.close()

for test in tests:
    print(test.shape)

for test in tests:
    plt.figure()
    plt.imshow(np.uint8(test[0]), interpolation='none')
    plt.show()
